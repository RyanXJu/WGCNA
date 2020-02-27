## 

boxplot(Overall_Survival_Time_days~Sex, data=datTraits, main="OS vs Sex",
        xlab="sex", ylab="OS days")


boxplot(Overall_Survival_Time_days~dx_FAB, data=datTraits_FAB, main="OS vs FAB",
        xlab="FAB", ylab="OS days")


boxplot(Overall_Survival_Time_days~WHO.2008, data=datTraits_WHO, main="OS vs WHO",
        xlab="WHO", ylab="OS days")
